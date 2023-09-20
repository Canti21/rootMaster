const burguerMenu = document.querySelector('.burguer-menu');
const mobileMenu = document.querySelector('.mobile-menu');
const darken = document.querySelector('.darken')

const inputFields = document.querySelectorAll('.input-field')
const errorMessages = document.querySelectorAll('.error')

burguerMenu.addEventListener('click', toggleMobileMenu);

function toggleMobileMenu()
{
    mobileMenu.classList.toggle('inactive')
    darken.classList.toggle('inactive')
}

inputFields.forEach(function (inputField, index)
{
    inputField.addEventListener("input", function ()
    {
        if (!inputField.checkValidity())
        {
            errorMessages[index].classList.remove("inactive")
            inputField.classList.add("error-input"); // Add the red border class
        }
        else
        {
            errorMessages[index].classList.add("inactive");
            inputField.classList.remove("error-input"); // Remove the red border class
        }
    });
});

function validateForm()
{
    let isValid = true;

    inputFields.forEach(function (inputField, index)
    {
        if (!inputField.checkValidity())
        {
            errorMessages[index].classList.remove("inactive")
            inputField.classList.add("error-input"); // Add the red border class
            isValid = false;
        }
        else
        {
            errorMessages[index].innerHTML = "";
        }
    });

    if (isValid)
    {
        // If all fields are valid, manually submit the form.
        document.querySelector("form").submit();
    }
}